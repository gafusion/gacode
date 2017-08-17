!---------------------------------------------------------------
! gyro_collision_main.f90 
!
! PURPOSE:
!  Collision step for Lorentz collisions.
!---------------------------------------------------------------

subroutine gyro_collision_main

  use gyro_globals
  use gyro_pointers
  use ompdata

  !----------------------------------------
  implicit none
  !
  integer :: ic
  integer :: n_i
  integer :: n_j
  integer :: n_k
  integer :: n_d
  integer, external :: parallel_dim
  !----------------------------------------

  !---------------------
  ! BEGIN Collision step
  !---------------------

  do ic=1,n_coll

     call gyro_timer_in('Coll.-step')

     is = c_map(ic)

     ! H_a = h_a + Psi_a
     !
     ! d(H_a)/dt = C(H_a) + d(Psi_a)/dt
     !
     ! We will ignore d(Psi_a)/dt

!$omp parallel private(p_nek_loc,ie,k,p_nek,m,i)
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        do m=1,n_stack
           do i=ibeg,iend

              f_coll(m,i,p_nek_loc) = h(m,i,p_nek_loc,is)+&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)

           enddo ! i
        enddo ! m

     enddo ! p_nek
!$omp end parallel
     !-----------------------------------------------------

     !-----------------------------------------------
     ! Step 3:
     ! 
     ! Implicit collision solve:
     !
     ! Permute p_nek,i -> p_ine,k
     n_i = n_n_1*n_energy
     n_j = n_lambda
     n_k = n_x
     n_d = parallel_dim(n_k*n_i,n_proc_1)
     !
     call gyro_timer_out('Coll.-step')
     call gyro_timer_in('Coll.-comm')
     !
     call rTRANSP_INIT(n_i,n_j,n_k,n_stack,NEW_COMM_1)
     call rTRANSP_DO(f_coll,h_c)
     call rTRANSP_CLEANUP
     !
     call gyro_timer_out('Coll.-comm')
     call gyro_timer_in('Coll.-step')
     !
     ! Crank-Nicholson advance of H.
     !
     call gyro_collision_kernel(ic)
     !
     ! Permute p_ine,k -> p_nek,i
     n_i = n_x
     n_j = n_n_1*n_energy
     n_k = n_lambda
     n_d = parallel_dim(n_j*n_k,n_proc_1)
     !
     call gyro_timer_out('Coll.-step')
     call gyro_timer_in('Coll.-comm')
     !
     call fTRANSP_INIT(n_i,n_j,n_k,n_stack,NEW_COMM_1)
     call fTRANSP_DO(h_c,fb_coll)
     call fTRANSP_CLEANUP
     !
     call gyro_timer_out('Coll.-comm')
     call gyro_timer_in('Coll.-step')
     !
     ! NOW:
     !
     ! f_coll  = H = h+Psi
     ! fb_coll = H_new 
     !------------------------------------------------

     !-------------------------------------------------------- 
     ! Step 4:
     !
     ! Correct for small amount of particle nonconservation.
     ! This avoids a weak long-wavelength instability for n=0.
     !
     if (is < indx_e .and. electron_method /= 3) then
        ! Total conservation
        fb_coll(:,:,:) = fb_coll(:,:,:)-f_coll(:,:,:)
        call gyro_conserve_all
        h(:,:,:,is) =  h(:,:,:,is)+fb_coll(:,:,:)
     endif
     if (is == indx_e) then 
        ! Number conservation
        h(:,:,:,is) = h(:,:,:,is)+fb_coll(:,:,:)-f_coll(:,:,:)
        call gyro_conserve_number
     endif
     !
     !--------------------------------------------------------

     !---------------------
     ! END Collision step
     !---------------------

     call gyro_timer_out('Coll.-step')

  enddo

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_collision_main done]'
  endif

end subroutine gyro_collision_main
