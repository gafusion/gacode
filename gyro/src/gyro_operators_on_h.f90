!-----------------------------------------------------
! gyro_operators_on_h.f90
!
! PURPOSE:
!  Given h, compute the gyroaverages for all gyrokinetic 
!  species, and limiting forms for DK electrons.
!
!  This is valid for periodic OR nonperiodic boundary
!  conditions. 
!
!  NOTE:
!   This routine is expensive so some coding is inelegant 
!   for the sake of speed.
!-----------------------------------------------------

subroutine gyro_operators_on_h

  use gyro_globals
  use gyro_pointers

  !-------------------------------------------------------
  implicit none
  !
  complex :: temp
  complex :: temp2
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: hh
  !-------------------------------------------------------

  call gyro_timer_in('Gyroave-h')

  !-------------------------------------------------------
  ! Gyrokinetic species:
  !
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do is=1,n_gk

        do i=1,n_x
           hh(:,i) = h(:,i,p_nek_loc,is)
        enddo
        if (boundary_method == 1) then
           do i=1-m_gyro,0
              hh(:,i) = h(:,i+n_x,p_nek_loc,is)
           enddo
           do i=n_x+1,n_x+m_gyro-i_gyro
              hh(:,i) = h(:,i-n_x,p_nek_loc,is)
           enddo
        else
           do i=1-m_gyro,0
              hh(:,i) = 0.0
           enddo
           do i=n_x+1,n_x+m_gyro
              hh(:,i) = 0.0
           enddo
        endif

        if (n_field < 3) then

!$omp parallel do private(m,i_diff,temp)
           do i=1,n_x
              do m=1,n_stack
                 temp = (0.0,0.0)
                 do i_diff=-m_gyro,m_gyro-i_gyro
                    temp = temp+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)*hh(m,i+i_diff)
                 enddo ! i_diff
                 gyro_h(m,i,p_nek_loc,is) = temp
              enddo ! m
           enddo ! i
!$omp end parallel do

        else

!$omp parallel do private(i_diff,m,temp,temp2)
           do i=1,n_x
              do m=1,n_stack
                 temp = (0.0,0.0)
                 temp2 = (0.0,0.0)
                 do i_diff=-m_gyro,m_gyro-i_gyro
                    temp = temp+&
                         w_gyro0(m,i_diff,i,p_nek_loc,is)*hh(m,i+i_diff)
                    temp2 = temp2+&
                         w_gyro1(m,i_diff,i,p_nek_loc,is)*hh(m,i+i_diff)
                 enddo ! i_diff
                 gyro_h(m,i,p_nek_loc,is) = temp
                 gyro_h_aperp(m,i,p_nek_loc,is) = temp2
              enddo ! m
           enddo ! i
!$omp end parallel do

        endif

     enddo ! is
  enddo ! p_nek
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Drift-kinetic electrons:
  !
  ! J0        -> 1
  ! (J0+J2)/2 -> 1/2
  !
  if (electron_method == 2) then
     gyro_h(:,:,:,n_spec) = h(:,:,:,n_spec)
     if (n_field > 2) gyro_h_aperp(:,:,:,n_spec) = 0.5*h(:,:,:,n_spec)
  endif
  !-------------------------------------------------------

  call gyro_timer_out('Gyroave-h')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_operators_on_h done]'
  endif

end subroutine gyro_operators_on_h
