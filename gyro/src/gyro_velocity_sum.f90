!--------------------------------------------------------------
! get_velocity_sum.f90
!
! PURPOSE:
!  Evaluation kernel for RHS velocity-space (source) integrals.
!--------------------------------------------------------------

subroutine gyro_velocity_sum(field)

  use mpi
  use gyro_globals
  use gyro_pointers
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: field
  complex, dimension(n_blend,n_x) :: sum_loc
  complex, dimension(n_blend,n_x) :: sum_glob
  complex, dimension(n_stack) :: gz
  !---------------------------------------------------

  call gyro_timer_in('Velocity-sum')

  !----------------------------------------------------
  ! Now, compute blending projections:
  !
!$omp parallel private(p_nek_loc,p_nek,k,ck,gz,m0,is,m,j,i)
  sum_loc(:,ibeg:iend) = (0.0,0.0)
  !
  select case (field)

  case (1)

     ! Phi
     !
     ! sum_s FV[(F*_j) z_s*<hi>]
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        k  = nek_k(p_nek)   
        ck = class(k)

        do i = ibeg, iend

           gz(:) = (0.0,0.0)
           do is=1,n_kinetic
              gz(:) = gz(:)+z(is)*gyro_h(:,i,p_nek_loc,is)
           enddo

           do m=1,n_stack
              m0 = m_phys(ck,m)
              do j=1,n_blend
                 sum_loc(j,i) = sum_loc(j,i)+gz(m)*&
                      cs_blend(j,m0,i,p_nek_loc)
              enddo
           enddo ! m
        enddo 

     enddo 

  case (2)

     ! A_parallel
     !
     ! sum_s FV[(F*_j) z_s*v_s*<h_s>]

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   
        ck = class(k)

        do i = ibeg, iend

           gz(:) = (0.0,0.0)
           do is=1,n_kinetic
              gz(:) = gz(:)+z(is)*gyro_h(:,i,p_nek_loc,is)*&
                   v_para(:,i,p_nek_loc,is)
           enddo

           do m=1,n_stack
              m0 = m_phys(ck,m)
              do j=1,n_blend
                 sum_loc(j,i) = sum_loc(j,i)+gz(m)*&
                      cs_blend(j,m0,i,p_nek_loc)
              enddo
           enddo ! m
        enddo 

     enddo

  case (3)

     ! B_parallel
     !
     ! sum_s FV[(F*_j) G_perp(hi)*(-T_s*ene*lambda)]

     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   
        ck = class(k)

        do i = ibeg, iend

           gz(:) = (0.0,0.0)
           do is=1,n_kinetic 
              gz(:) = gz(:)-gyro_h_aperp(:,i,p_nek_loc,is)*&
                   tem_s(is,i)*energy(ie)*lambda(i,k)
           enddo

           do m=1,n_stack
              m0 = m_phys(ck,m)
              do j=1,n_blend
                 sum_loc(j,i) = sum_loc(j,i)+gz(m)*&
                      cs_blend(j,m0,i,p_nek_loc)
              enddo
           enddo ! m

        enddo 

     enddo 

  end select
!$omp end parallel 

  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Final velocity-space sum on subgroup
  !
  call MPI_ALLREDUCE(sum_loc,&
       sum_glob,&
       size(sum_glob),&
       MPI_DOUBLE_COMPLEX,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)
  !--------------------------------------------------------------

  if (n_1(in_1) == 0 .and. boundary_method == 1) then

     ! Eliminate average in periodic case
     sum_glob(:,n_x) = (0.0,0.0)

  endif

  !--------------------------------------------------------------
  select case (field)

  case (1) 
     vel_sum_p = sum_glob
  case (2) 
     vel_sum_a = sum_glob
  case (3) 
     vel_sum_aperp = sum_glob

  end select
  !--------------------------------------------------------------

  call gyro_timer_out('Velocity-sum')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_velocity_sum done]'
  endif

end subroutine gyro_velocity_sum
