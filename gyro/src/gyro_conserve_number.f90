!------------------------------------------------
! gyro_conserve_number.f90
!
! PURPOSE:
!  Enforce number conservation by correcting h 
!  (based on value in h_old) to remove the 
!  component which alters number
!------------------------------------------------

subroutine gyro_conserve_number

  use mpi
  use gyro_globals
  use gyro_pointers
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_blend,n_x) :: vel_sum_loc
  complex, dimension(n_blend,n_x) :: vel_sum
  !---------------------------------------------------

  !----------------------------------------------------
  ! Now, compute blending projections:
  !
  !
!$omp parallel private(p_nek_loc,ie,k,ck,m0,p_nek,i,m)
  vel_sum_loc(:,ibeg:iend) = (0.0,0.0)

  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=ibeg,iend

        do m=1,n_stack

           m0 = m_phys(ck,m)

           ! Poisson Equation (RHS)
           !
           ! FV[(F*_j) he]

           vel_sum_loc(:,i) = vel_sum_loc(:,i)+ &
                (h_old(m,i,p_nek_loc,is)-h(m,i,p_nek_loc,is))* &
                cs_blend(:,m0,i,p_nek_loc)

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
!$omp end parallel
  !--------------------------------------------------------------

  call MPI_ALLREDUCE(vel_sum_loc,&
       vel_sum,&
       n_blend*n_x,&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  !---------------------------------------------------
  ! Solve for blending coefficients
  !
!$omp parallel do private(info) schedule(static)
  do i=1,n_x

     call ZGETRS('N',&
          n_blend,&
          1,&
          ff_mm(:,:,i,2),&
          n_blend,&
          ff_mm_piv(:,i,2),&
          vel_sum(:,i),&
          n_blend,&
          info)

  enddo ! i
  !---------------------------------------------------


  !---------------------------------------------------
  ! Correct distribution
  !
!$omp parallel private(p_nek_loc,ie,k,ck,m0,is,p_nek,i,m)
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=ibeg,iend

        do m=1,n_stack

           m0 = m_phys(ck,m)

           is = n_kinetic

           h(m,i,p_nek_loc,is) = h(m,i,p_nek_loc,is)+ &
                sum(vel_sum(:,i)*c_blend(:,m0,i,p_nek_loc))

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc
!$omp end parallel

end subroutine gyro_conserve_number
