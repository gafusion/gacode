!-----------------------------------------------------
! get_delta_he.f90
!
! PURPOSE:
!  Compute the "explicit" part of h_e.
!-----------------------------------------------------

subroutine gyro_get_delta_he

  use gyro_globals
  use gyro_pointers
  use ompdata

  !------------------------------------------
  implicit none
  !
  integer :: n_s2
  integer :: mp
  !
  complex :: h_temp(n_stack,n_x,n_nek_loc_1)
  !------------------------------------------

  n_s2 = n_stack/2

!$omp parallel private(p_nek_loc,k,ck,p_nek,i,mp,m)
  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     k  = nek_k(p_nek)   

     ck = class(k)

     h_temp(:,ibeg:iend,p_nek_loc)   = h(:,ibeg:iend,p_nek_loc,n_spec)
     h(:,ibeg:iend,p_nek_loc,n_spec) = (0.0,0.0)

     if (ck == 1) then

        do i = ibeg, iend
           do mp=1,n_s2
              do m=1,n_s2
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
           do mp=n_s2+1,n_stack
              do m=n_s2+1,n_stack
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
        enddo ! i

     else

        do i = ibeg, iend
           do mp=1,n_stack
              do m=1,n_stack
                 h(m,i,p_nek_loc,n_spec) = h(m,i,p_nek_loc,n_spec)+&
                      o_advect(m,mp,i,p_nek_loc)*h_temp(mp,i,p_nek_loc)
              enddo ! m
           enddo ! mp
        enddo ! i

     endif

  enddo ! p_nek_loc
!$omp end parallel


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_get_delta_he done]'
  endif

end subroutine gyro_get_delta_he
